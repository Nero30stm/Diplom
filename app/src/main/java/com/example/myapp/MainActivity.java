package com.example.myapp;

import android.os.Environment;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.app.Notification;
import android.app.NotificationManager;
import android.app.PendingIntent;
import java.util.List;
import android.app.ListActivity;
import android.app.ProgressDialog;
import android.content.ActivityNotFoundException;
import android.content.Intent;
import android.content.pm.ApplicationInfo;
import android.content.pm.PackageManager;
import android.os.AsyncTask;
import android.view.View;
import android.widget.ListView;
import android.widget.Toast;

import java.util.ArrayList;
import android.media.MediaRecorder;
import android.content.Context;
import android.content.res.Resources;
import android.graphics.BitmapFactory;
import android.support.v4.app.NotificationCompat;
import android.app.Service;
import android.app.Activity;
import java.io.File;
import android.media.MediaPlayer;

public class MainActivity extends ListActivity {

    //private ListView appList;

    private PackageManager packageManager = null;
    private List <ApplicationInfo> applist = null;
    private AppListAdapter listadapter = null;

    private MediaRecorder mediaRecorder;
    private MediaPlayer mediaPlayer;
    private String fileName;
    // Идентификатор уведомления
    private static final int NOTIFY_ID = 101;
    private NotificationManager nm;
    MediaRecorder myAudioRecorder = new MediaRecorder();

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

       // appList = (ListView) findViewById(R.id.appList);

        //PackageManager pm = getPackageManager();
        //List<ApplicationInfo> packages = pm
          //      .getInstalledApplications(PackageManager.GET_META_DATA);

        //AppListAdapter appListAdapter = new AppListAdapter(this, pm, packages);
        //appList.setAdapter(appListAdapter);





        packageManager = getPackageManager();

        new LoadApplications().execute();

        nm = (NotificationManager) getApplicationContext().getSystemService(Context.NOTIFICATION_SERVICE);
        showNotification();

    }




    @Override
    protected void onListItemClick(ListView l, View v, int position, long id) {
        super.onListItemClick(l, v, position, id);

        ApplicationInfo app = applist.get(position);


        Intent intent = new Intent(this, SpeechRecognition.class);
        intent.putExtra("packageName", app.packageName);
        startActivity(intent);
       /* try{
            Intent intent = packageManager.getLaunchIntentForPackage(app.packageName);

            if(intent != null) {
                startActivity(intent);
            }
        } catch(ActivityNotFoundException e) {
            Toast.makeText(MainActivity.this, e.getMessage(), Toast.LENGTH_LONG).show();
        } catch(Exception e) {
            Toast.makeText(MainActivity.this, e.getMessage(), Toast.LENGTH_LONG).show();
        }*/

    }

    private List <ApplicationInfo> checkForLaunchIntent(List <ApplicationInfo> list) {

        ArrayList <ApplicationInfo> appList = new ArrayList <ApplicationInfo> ();

        for(ApplicationInfo info : list) {
            try{
                if(packageManager.getLaunchIntentForPackage(info.packageName) != null) {
                    appList.add(info);
                }
            } catch(Exception e) {
                e.printStackTrace();
            }
        }

        return appList;
    }

    private class LoadApplications extends AsyncTask<Void, Void, Void> {

        private ProgressDialog progress = null;

        @Override
        protected Void doInBackground(Void... params) {

            applist = checkForLaunchIntent(packageManager.getInstalledApplications(PackageManager.GET_META_DATA));

            listadapter = new AppListAdapter(MainActivity.this, R.layout.app_item, applist);

            return null;
        }

        @Override
        protected void onPostExecute(Void result) {
            setListAdapter(listadapter);
            progress.dismiss();
            super.onPostExecute(result);
        }

        @Override
        protected void onPreExecute() {
            progress = ProgressDialog.show(MainActivity.this, null, "Loading apps info...");
            super.onPreExecute();
        }
    }


    public void showNotification () {
        Notification.Builder builder = new Notification.Builder(getApplicationContext());
        //Intent notificationIntent = getPackageManager().getLaunchIntentForPackage("com.google.android.youtube");
        Intent notificationIntent = new Intent(this, SpeechRecognition.class);
        PendingIntent contentIntent = PendingIntent.getActivity(this,0, notificationIntent, PendingIntent.FLAG_CANCEL_CURRENT);
        builder.setContentIntent(contentIntent)
                // обязательные настройки
                .setSmallIcon(R.drawable.ic_mic)
                //.setContentTitle(res.getString(R.string.notifytitle)) // Заголовок уведомления
                .setContentTitle("Голосовая команда")
                //.setContentText(res.getString(R.string.notifytext))
                .setContentText("Нажмите для распознавания команды") // Текст уведомления
                // необязательные настройки
                .setLargeIcon(BitmapFactory.decodeResource(getApplicationContext().getResources(), R.drawable.ic_micro)) // большая
                // картинка
                //.setTicker(res.getString(R.string.warning)) // текст в строке состояния
                .setTicker("По")
                .setWhen(System.currentTimeMillis());
                //.setAutoCancel(true); // автоматически закрыть уведомление после нажатия

        NotificationManager notificationManager = (NotificationManager) getSystemService(Context.NOTIFICATION_SERVICE);
        // Альтернативный вариант
        // NotificationManagerCompat notificationManager = NotificationManagerCompat.from(this);
        nm.notify(NOTIFY_ID, builder.build());
    }
}


/*public class MainActivity extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
    }
}*/

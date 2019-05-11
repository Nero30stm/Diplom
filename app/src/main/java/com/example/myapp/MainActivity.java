package com.example.myapp;

import android.os.Environment;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.app.Notification;
import android.app.NotificationManager;
import android.app.PendingIntent;
import java.util.List;
import android.media.MediaRecorder;
import android.app.Notification;
import android.app.NotificationManager;
import android.app.PendingIntent;
import android.content.Context;
import android.content.Intent;
import android.content.res.Resources;
import android.graphics.BitmapFactory;
import android.os.Bundle;
import android.support.v4.app.NotificationCompat;
import android.support.v7.app.AppCompatActivity;
import android.view.View;
import android.content.Intent;
import android.app.Service;
import android.app.Activity;
import android.content.pm.ApplicationInfo;
import android.content.pm.PackageManager;
import java.io.File;
import android.view.View;
import android.widget.ListView;
import android.media.MediaPlayer;
import android.widget.Toast;
import android.app.Activity;
import android.media.MediaPlayer;
import android.media.MediaRecorder;
import android.os.Bundle;
import android.os.Environment;

public class MainActivity extends Activity {

    private ListView appList;
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

        appList = (ListView) findViewById(R.id.appList);

        PackageManager pm = getPackageManager();
        List<ApplicationInfo> packages = pm
                .getInstalledApplications(PackageManager.GET_META_DATA);

        AppListAdapter appListAdapter = new AppListAdapter(this, pm, packages);
        appList.setAdapter(appListAdapter);
        

        nm = (NotificationManager) getApplicationContext().getSystemService(Context.NOTIFICATION_SERVICE);
        showNotification();

    }

    public void onMyButtonClick(View view)
    {
        //intent.putExtra("name", name.getText().toString());

        /*Intent sendIntent = getPackageManager().getLaunchIntentForPackage("com.google.android.youtube");

        if (sendIntent != null) {
            sendIntent.addFlags(Intent.FLAG_ACTIVITY_NEW_TASK);
            startActivity(sendIntent);
        }*/

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
